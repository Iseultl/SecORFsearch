process CONCATENATE_OG {
    input:
    path og_files

    output:
    path "concatenated_og.txt"

    script:
    """
    #!/bin/bash
    echo ${og_files}
    cat ${og_files.join(' ')} > concatenated_og.txt
    """
}