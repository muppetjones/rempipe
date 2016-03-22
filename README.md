

## CRITICAL

* [Known Python 3.4-6 bug](https://bugs.python.org/issue18622)

    This bug is triggered in `tests.decorators.test_io`. The offending test
    is skipped by default, but may be reenabled after patching Python. The
    required patch is included in the `patch` directory.

    To apply the patch, run the following command:
    ```shell
    patch -p2 -b < file.patch
    ```
    where `file.patch` is the patch file.
