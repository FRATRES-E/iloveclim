#!/usr/bin/env bash
# =============================================================================
# Hooks de la composante CLIO
# =============================================================================
# Vit avec le code de couplage CLIO. Le moteur le source si le descripteur
# déclare un hook.
#
# clio_tracers : dérive le nombre de traceurs océan (nsmax) et copie la version
# correspondante de netcdfout.param. Remplace le système de sentinelles
# (.activated) et le choix manuel netcdfout.param-scalNN.
#
# SOURCE DE VÉRITÉ : para0_mod.f90 calcule nsmax par préprocesseur à partir de
# choixcomposantes.h. On ne duplique PAS cette logique ; on la LIT. Le nom du
# fichier de sortie se calcule (netcdfout.param-scal${nsmax}), il ne se choisit
# pas par une table de cas.
#
# PRIORITÉ -F (local à ce hook) : si l'utilisateur a fourni via -F un para0_mod.f90
# ou un choixcomposantes.h modifié, CE SONT EUX qui font autorité (cohérence avec
# le binaire qui sera compilé). Le hook applique cette priorité pour ses 2 fichiers
# sans toucher à la logique de copie globale du script.
#
# Contrat : reçoit
#   $1 = run_dir
#   $2 = target (répertoire où déposer netcdfout.param, relatif à run_dir ; ""=run_dir)
# Utilise les variables globales du script : comp_dir, emic_dir, source_dir,
# OCEAN, verbose, run_label, cpp (préprocesseur).
# Retourne 0 si succès, non-zéro si dérivation impossible (arrêt du run voulu).
# =============================================================================

clio_tracers() {
    local run_dir=$1
    local target_sub=$2
    local dest="${run_dir}/${target_sub}"
    [ -z "$target_sub" ] && dest="${run_dir}"

    # --- Préprocesseur : réutiliser celui du Makefile si fourni, sinon cpp ----
    local PP="${CPP_FOR_TRACERS:-cpp}"
    if ! command -v "$PP" >/dev/null 2>&1; then
        echo "[clio_tracers] ERROR: préprocesseur '$PP' introuvable." >&2
        return 1
    fi

    # --- Fichier à préprocesser --------------------------------------------------
    # On préprocesse SUR PLACE dans comp_dir/sources : c'est là que se trouvent,
    # à ce stade, les versions qui FONT AUTORITÉ. En effet la surcharge -F (si
    # présente) a DÉJÀ déposé ses para0_mod.f90 / choixcomposantes.h dans
    # comp_dir/sources en écrasant les versions de base. Préprocesser ici :
    #   (1) reproduit exactement l'environnement d'include de la compilation
    #       (tous les .h voisins présents, quel que soit le style de #include) ;
    #   (2) garantit la cohérence nsmax <-> binaire réellement compilé ;
    #   (3) évite l'artefact d'un tmp isolé où un #include indirect échouait.
    local src_dir="${comp_dir}/sources"
    local para0="${src_dir}/para0_mod.f90"

    if [ ! -f "$para0" ]; then
        echo "[clio_tracers] ERROR: para0_mod.f90 introuvable ($para0)." >&2
        return 1
    fi

    # --- Extraire nsmax (robuste : trailing ws, commentaires !, nsmax_TS) -----
    # Préprocessing depuis src_dir (-I sur src_dir) pour résoudre les #include
    # comme le fera le compilateur.
    local nsmax
    nsmax=$( cd "$src_dir" && "$PP" -P -traditional-cpp -I. para0_mod.f90 2>/dev/null \
        | grep -E 'parameter[[:space:]]*::[[:space:]]*nsmax[[:space:]]*=' \
        | grep -vE 'nsmax_' \
        | head -1 \
        | sed -E 's/.*nsmax[[:space:]]*=[[:space:]]*([0-9]+).*/\1/' )

    if ! [[ "$nsmax" =~ ^[0-9]+$ ]]; then
        echo "[clio_tracers] ERROR: extraction de nsmax échouée (obtenu: '${nsmax}')." >&2
        echo "[clio_tracers]        fichier: $para0" >&2
        echo "[clio_tracers]        (vérifier le #include de choixcomposantes.h et le préprocesseur)" >&2
        return 1
    fi

    # --- Copier la version de netcdfout.param correspondante ------------------
    local src_param="${emic_dir}/param/${OCEAN}/netcdfout.param-scal${nsmax}"
    if [ ! -f "$src_param" ]; then
        echo "[clio_tracers] ERROR: nsmax=${nsmax} dérivé, mais fichier absent :" >&2
        echo "[clio_tracers]        $src_param" >&2
        echo "[clio_tracers]        (créer netcdfout.param-scal${nsmax} ou vérifier la config)" >&2
        return 1
    fi

    cp "${src_param}" "${dest}/netcdfout.param" || {
        echo "[clio_tracers] ERROR: copie de netcdfout.param échouée vers ${dest}." >&2
        return 1
    }

    [ "${verbose:-0}" -ge 1 ] && \
        echo "  X CLIO: ${nsmax} traceurs océan -> netcdfout.param-scal${nsmax}"
    return 0
}
