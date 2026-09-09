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
# injecte que ce qui varie : num_years, start_year, NSKIP. Séparation nette
# entre le "quoi" (template) et le "comment" (ce hook).
#
# Contrat : reçoit
#   $1 = run_dir
#   $2 = target (nom du namelist à produire, relatif à run_dir ; ""=namelist)
# Variables runtime utilisées (globales du script) : num_years, start_year, NSKIP.
# Le template est cherché à côté du cache (copié là par le générateur) ; à défaut
# dans le répertoire source de la composante.
# Retourne 0 si succès, non-zéro sinon (template ou variable manquant => arrêt).
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

    # --- Substituer les marqueurs et écrire le namelist -----------------------
    # sed sur chaque marqueur. Les valeurs sont numériques (pas de caractère
    # spécial sed attendu), mais on reste prudent sur le délimiteur.
    sed -e "s|@NUM_YEARS@|${num_years}|g" \
        -e "s|@START_YEAR@|${start_year}|g" \
        -e "s|@NSKIP@|${NSKIP}|g" \
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
