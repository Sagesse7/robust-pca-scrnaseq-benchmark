# Computational environments

The formal large-scale analyses ran under R 4.4.1 on x86_64 Rocky Linux. The
corresponding records are:

- `formal_server_sessionInfo.txt`;
- `formal_server_packages.csv`;
- `formal_server_conda_environment.yml`.

Figure generation and local package verification used R 4.5.0 on macOS. Those
records are:

- `figure_generation_sessionInfo.txt`;
- `figure_generation_package_versions.csv`.

CPU-based execution was used throughout. Numerical-library threading was fixed
to one thread per method process for the formal simulation and benchmark jobs;
real-data runners used the explicitly recorded worker count. No GPU is required.
