

Channel.fromPath("tests/*", type: "any").set { to_test }

process run_tests {
    container "$params.annotation_container"

    input:
        path tests_dir from to_test.collect()

    script:
    """
        ./test.py
    """
}
