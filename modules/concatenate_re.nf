process CONCATENATE_RE {
    input:
    path re_files

    output:
    path "concatenated_re.txt"

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail
    echo ${re_files}
    cat ${re_files.join(' ')} > concatenated_re.txt
    """
}