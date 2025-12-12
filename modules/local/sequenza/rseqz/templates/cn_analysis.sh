#!/usr/bin/env bash
set -euo pipefail

# Compute variables inside template
prefix="${meta.id}_${purity}"
seqz=\$(echo "${seqz_bin}" | sed 's/\\.gz\$//')
if [ ! -d \${prefix} ]; then
  echo "mkdir \${prefix}"
  mkdir \${prefix}
fi

# zcat ${seqz_bin} > \${seqz}
Rscript ${projectDir}/templates/analyse_cn_sequenza.Rscript \
${seqz_bin} \${prefix} ${meta.id} ${meta.sex} ${meta.ploidy} ${purity} ${meta.gamma}

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    sequenza: 2.1.1
END_VERSIONS
