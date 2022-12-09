#!/usr/bin/env python3

import os
import shutil

import annotations

class TestCheckGbk:
    def __init__(self, test_name, test_dir, run):
        self.name = test_name
        self.test_dir = test_dir
        self.passed = False
        self.run = run

    def cleanup(self):
        shutil.remove(self.test_dir)


def check_gbk_run_test_1(test):
    annotations.check_gbk(test.test_dir + "/test.csv")

def check_gbk_run_test_1(test):
    annotations.check_gbk(test.test_dir + "/test.csv")


tests_list = [
    TestCheckGbk("Locus tag collision", "test1", run=check_gbk_run_test_1),
    TestCheckGbk("Organism collision", "test2", run=check_gbk_run_test_2),
    TestCheckGbk("Different organism in same gbk", "test3"),
    TestCheckGbk("Organism renaming", "test4"),
    TestCheckGbk("Organism renaming, no name column", "test5"),
    TestCheckGbk("Organism renaming, no names", "test6"),
    TestCheckGbk("Wrong format", "test7")
]


failed = []
n_tests = 0
for test in tests_list:
    test.run(test)

    if not test.passed:
        failed.append(test)
    n_tests += 1

if len(failed) == 0:
    print(f"All tests passed ({n_tests}/{n_tests})")
else:
    print(f"{len(failed)}/{n_tests} tests have failed:")
    for failed_test in failed:
        print(failed_test.name)

print("Cleaning up")
for test in tests_list:
    test.cleanup()
