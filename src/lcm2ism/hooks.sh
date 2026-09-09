#!/usr/bin/env bash
# =============================================================================
# Hooks de la composante GRISLI
# =============================================================================
# Vit avec le code de la composante (src/grisli/). Le générateur le recopie
# à côté du cache ; le moteur le source si le descripteur déclare un hook.
#
# Contrat d'un hook restart : reçoit
#   $1 = run_dir   (racine du run)
#   $2 = target    (fichier de paramètres à patcher, relatif à run_dir)
# Retourne 0 en cas de succès. La logique est libre : c'est du bash ordinaire.
# =============================================================================

grisli_restart() {
    local run_dir=$1
    local target="${run_dir}/$2"
    local fname

    if [ ! -f "$target" ]; then
        echo "[grisli_restart] ERROR: target not found: $target" >&2
        return 1
    fi

    if compgen -G "${run_dir}/startdata/*grestart*.nc" >/dev/null 2>&1; then
        fname=$(basename $(ls -1 ${run_dir}/startdata/*grestart*.nc | tail -n 1))
        sed -i "s%THERESTARTFILEGRISLI%./startdata/${fname}%" "$target"
        sed -i "s/THERESTARTCOMPTEUR/1/;s/THERESTARTIOUT/2/" "$target"
        echo "  X GRISLI will use a nc restart file (ALL variables): ${fname}"
    elif compgen -G "${run_dir}/startdata/*grestart*.cptr" >/dev/null 2>&1; then
        fname=$(basename $(ls -1 ${run_dir}/startdata/*grestart*.cptr | tail -n 1))
        sed -i "s%THERESTARTFILEGRISLI%./startdata/${fname}%" "$target"
        sed -i "s/THERESTARTCOMPTEUR/1/;s/THERESTARTIOUT/2/" "$target"
        echo "  X GRISLI will use a cptr restart file (ALL variables): ${fname}"
    else
        sed -i "s/THERESTARTCOMPTEUR/0/;s/THERESTARTIOUT/2/" "$target"
        echo "  X GRISLI: no restart file provided, starting from scratch"
    fi
    return 0
}
