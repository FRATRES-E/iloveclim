#!/usr/bin/env bash
# =============================================================================
# Hooks de la composante CORE (socle commun)
# =============================================================================
# Vit avec le descripteur core (src/core/). Le moteur le source si le
# descripteur déclare un hook.
# =============================================================================

# core_namelist : génère le namelist principal du run à partir d'un template
# externe (namelist.template), en substituant les trois variables du run.
#
# Le template est un fichier INERTE et lisible (contenu stable) ; le hook n'y
# injecte que ce qui varie : num_years, start_year, NSKIP, et les métadonnées
# netCDF du groupe ncmeta (institution, author, source). Séparation nette
# entre le "quoi" (template) et le "comment" (ce hook).
#
# Métadonnées netCDF (groupe ncmeta, global attributes CF des fichiers io_nc) :
#   @NC_INSTITUTION@ <- ${ILOVECLIM_INSTITUTION:-NotSet}  (réglé à l'installation)
#   @NC_AUTHOR@      <- ${USER:-${LOGNAME:-unknown}}
#   @NC_SOURCE@      <- "iLOVECLIM " + git describe --tags --always --dirty sur
#                       emic_dir ("unknown" si emic_dir n'est pas un dépôt git)
#
# Contrat : reçoit
#   $1 = run_dir
#   $2 = target (nom du namelist à produire, relatif à run_dir ; ""=namelist)
# Variables runtime utilisées (globales du script) : num_years, start_year, NSKIP.
# Le template est cherché à côté du cache (copié là par le générateur) ; à défaut
# dans le répertoire source de la composante.
# Retourne 0 si succès, non-zéro sinon (template ou variable manquant => arrêt).
# _core_nml_str : protège une chaîne libre pour l'insérer entre apostrophes
# dans un namelist Fortran ('  -> '') puis comme remplacement sed (\ & |).
_core_nml_str() {
    local v=$1
    v=${v//\'/\'\'}
    printf '%s' "$v" | sed -e 's/[\\&|]/\\&/g'
}

core_namelist() {
    local run_dir=$1
    local target_sub=$2
    [ -z "$target_sub" ] && target_sub="namelist"
    local dest="${run_dir}/${target_sub}"

    # --- Localiser le template ------------------------------------------------
    # Le générateur copie les fichiers auxiliaires du descripteur à côté du cache.
    # cache_dir est une globale du script ; on cherche d'abord là, puis dans le
    # répertoire source de la composante (mode dev/local).
    local tmpl=""
    if [ -n "${cache_dir:-}" ] && [ -f "${cache_dir}/core.namelist.template" ]; then
        tmpl="${cache_dir}/core.namelist.template"
    elif [ -f "${emic_dir}/defs/core/namelist.template" ]; then
        tmpl="${emic_dir}/defs/core/namelist.template"
    fi
    if [ -z "$tmpl" ]; then
        echo "[core_namelist] ERROR: template namelist introuvable" >&2
        echo "[core_namelist]        (ni ${cache_dir:-<cache_dir>}/core.namelist.template" >&2
        echo "[core_namelist]         ni ${emic_dir}/defs/core/namelist.template)" >&2
        return 1
    fi

    # --- Vérifier que les variables runtime sont définies ---------------------
    # Sémantique stricte : un namelist avec un marqueur non substitué (variable
    # vide) casserait le run silencieusement. On refuse plutôt.
    local v
    for v in num_years start_year NSKIP; do
        if [ -z "${!v:-}" ]; then
            echo "[core_namelist] ERROR: variable runtime '${v}' non définie ou vide" >&2
            return 1
        fi
    done

    # --- Métadonnées netCDF (texte libre) ---------------------------------------
    local nc_version nc_institution nc_author nc_source

        # --- Model version for the netCDF attribute "source" ----------------------
    local nc_vtag nc_hash nc_dirty
    nc_vtag=$(git -C "${emic_dir:-.}" log --format=%s 2>/dev/null | grep -m1 -oE '^v[0-9]+\.[0-9]+\.[0-9]+')
    if [ -n "$nc_vtag" ]; then
        nc_hash=$(git -C "${emic_dir:-.}" rev-parse --short HEAD 2>/dev/null)
        git -C "${emic_dir:-.}" diff --quiet HEAD 2>/dev/null || nc_dirty=", dirty"
        nc_version="${nc_vtag} (git ${nc_hash}${nc_dirty})"
    else
        nc_version=$(git -C "${emic_dir:-.}" describe --tags --always --dirty 2>/dev/null) || nc_version="unknown"
        [ -z "$nc_version" ] && nc_version="unknown"
    fi

    [ -z "$nc_version" ] && nc_version="unknown"
    nc_institution=$(_core_nml_str "${ILOVECLIM_INSTITUTION:-NotSet}")
    nc_author=$(_core_nml_str "${USER:-${LOGNAME:-unknown}}")
    nc_source=$(_core_nml_str "iLOVECLIM ${nc_version}")
    # --- Substituer les marqueurs et écrire le namelist -----------------------
    # sed sur chaque marqueur. Les valeurs numériques n'ont pas de caractère
    # spécial ; les chaînes libres sont protégées par _core_nml_str.
    sed -e "s|@NUM_YEARS@|${num_years}|g" \
        -e "s|@START_YEAR@|${start_year}|g" \
        -e "s|@NSKIP@|${NSKIP}|g" \
        -e "s|@NC_INSTITUTION@|${nc_institution}|g" \
        -e "s|@NC_AUTHOR@|${nc_author}|g" \
        -e "s|@NC_SOURCE@|${nc_source}|g" \
        "$tmpl" > "$dest" || {
        echo "[core_namelist] ERROR: écriture du namelist échouée ($dest)" >&2
        return 1
    }

    # --- Garde-fou : aucun marqueur ne doit subsister -------------------------
    if grep -q '@[A-Z_]*@' "$dest"; then
        echo "[core_namelist] ERROR: marqueur(s) non substitué(s) dans $dest :" >&2
        grep -o '@[A-Z_]*@' "$dest" | sort -u | sed 's/^/[core_namelist]   /' >&2
        return 1
    fi

    [ "${verbose:-0}" -ge 1 ] && \
        echo "  X CORE: namelist généré (nyears=${num_years}, irunlabel=${start_year}, nwrskip=${NSKIP})"
    return 0
}
