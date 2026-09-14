This is the initial submission of the package rsaber.

# Test environments

- Windows 11, R 4.5.0
- Windows 11, R 4.5.1
- Red Hat Enterprise Linux 8.10, R 4.3.0


# Results of: R CMD check --as-cran

0 errors | 0 warnings | 3 notes

There is one local NOTE on Windows:
* checking for future file timestamps ... unable to verify current time

This appears to be due to local system-clock verification being unavailable.
No files have future timestamps.

