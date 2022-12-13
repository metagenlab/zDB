#!/usr/bin/env python3

import os
import sys

import annotations

from ko_parsing import ko_tests


class TestCheckGbk:
    def __init__(self, test_name, test_dir, run):
        self.name = test_name
        self.test_dir = test_dir
        self.passed = False
        self.my_run = run

    def run(self):
        cur_dir = os.getcwd()
        os.chdir(self.test_dir)
        os.mkdir("filtered")
        self.my_run()
        os.chdir(cur_dir)


def check_gbk_run_test_1():
    annotations.check_gbk("test.csv")

def check_gbk_run_test_2():
    annotations.check_gbk("test.csv")

def check_gbk_run_test_3(test):
    annotations.check_gbk("test.csv")


tests_list = [
    TestCheckGbk("Locus tag collision different file", "test1", run=check_gbk_run_test_1),
    #TestCheckGbk("Organism collision", "test2", run=check_gbk_run_test_2),
    #TestCheckGbk("Different organism in same gbk", "test3"),
    #TestCheckGbk("Organism renaming", "test4"),
    #TestCheckGbk("Organism renaming, no name column", "test5"),
    #TestCheckGbk("Organism renaming, no names", "test6"),
    #TestCheckGbk("Locus tag collision same file", "test7", run=check_gbk_run_test_3),
    #TestCheckGbk("Wrong format", "test7"),
    *ko_tests.ko_tests
]


failed = []
n_tests = 0
for test in tests_list:
    test.run()

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

if len(failed)==0:
    sys.exit(0)
else:
    sys.exit("Some tests failed!")
