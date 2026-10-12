# syntax=docker/dockerfile:1
FROM debian:trixie-slim AS build

RUN apt-get update --error-on=any \
    && apt-get install -y --no-install-recommends \
        ca-certificates g++ git python3 zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Keep the compiler layer reusable when adding build tooling.
RUN apt-get update --error-on=any \
    && apt-get install -y --no-install-recommends cmake make \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY main.cpp CMakeLists.txt ./
COPY Kseq_Cpp/ ./Kseq_Cpp/

COPY tests/test_cli.py ./tests/test_cli.py

# CMake owns dependency download, compilation, tests and installation.
RUN --mount=type=secret,id=build_ca \
    --mount=type=cache,id=iecontig-seqan3,target=/build/_deps \
    if [ -f /run/secrets/build_ca ]; then export GIT_SSL_CAINFO=/run/secrets/build_ca; fi; \
    cmake -S . -B /build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON \
    && cmake --build /build --parallel 2 \
    && ctest --test-dir /build --output-on-failure -V \
    && cmake --install /build --prefix /out

FROM debian:trixie-slim AS runtime
RUN apt-get update --error-on=any \
    && apt-get install -y --no-install-recommends \
        ca-certificates libgomp1 libstdc++6 zlib1g \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /data \
    && chown 10001:10001 /data

COPY --from=build /out/ /usr/local/
ENV HOME=/tmp
WORKDIR /data
USER 10001:10001
ENTRYPOINT ["/usr/local/bin/IEContig"]
CMD ["--help"]
