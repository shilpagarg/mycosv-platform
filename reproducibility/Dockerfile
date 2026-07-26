# MycoSV core caller + Python analysis environment.
# Comparator tools and large reference databases (MycoMobilome, starbase) are
# fetched separately; see reproducibility/README.md.
FROM ubuntu:22.04

RUN apt-get update && apt-get install -y --no-install-recommends \
        g++ make python3 python3-pip samtools \
        ca-certificates curl \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir matplotlib numpy

WORKDIR /opt/mycosv
COPY main.cpp *.hpp ./
# Build the self-contained caller (no third-party tool in the calling path).
RUN g++ -O2 -std=c++17 -pthread main.cpp -o /usr/local/bin/fungi_graphsv_tol_bin

COPY reproducibility/ ./reproducibility/
# Build the element re-classifier used for the biology benchmark.
RUN g++ -O2 -std=c++17 -pthread -I. reproducibility/scripts/reannotate_elements.cpp \
        -o /usr/local/bin/reannotate_elements || true

ENTRYPOINT ["fungi_graphsv_tol_bin"]
CMD ["--help"]
