# On M1 Mac:
docker run -v "$(pwd)":"/opt/$(basename $(pwd))" --platform linux/amd64 -it --cap-add=SYS_PTRACE rocker/r-devel-san /bin/bash

# On Intel Mac:
docker run -v "$(pwd)":"/opt/$(basename $(pwd))" -it --cap-add=SYS_PTRACE rocker/r-devel-san /bin/bash

# To get httr package:
apt-get install libssl-dev
