

Channel.fromPath("tests/*").set { to_test }

process run_tests {
    container "$params.annotation_container"

    input:
        path tests_dir from to_test

    script:
    """
        ./test.py
    """
}
