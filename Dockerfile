# syntax=docker/dockerfile:1
FROM debian:trixie-slim AS build

RUN apt-get update --error-on=any \
    && apt-get install -y --no-install-recommends \
        ca-certificates g++ git python3 zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Pin SeqAn3 and use its recorded submodule commits.
RUN --mount=type=secret,id=build_ca \
    if [ -f /run/secrets/build_ca ]; then export GIT_SSL_CAINFO=/run/secrets/build_ca; fi; \
    git clone --depth 1 --branch 3.3.0 --recurse-submodules --shallow-submodules \
        https://github.com/seqan/seqan3.git /opt/seqan3 \
    && test "$(git -C /opt/seqan3 rev-parse HEAD)" = \
        6b681fb2eae5ab2997d293e99fc6a7f869a20316

WORKDIR /src
COPY main.cpp ./
COPY Kseq_Cpp/ ./Kseq_Cpp/

# Build only IEContig; the repository CMake file includes an absent example.
RUN mkdir -p /out \
    && g++ -std=c++20 -O2 -fopenmp -pthread \
        -I/opt/seqan3/include \
        -I/opt/seqan3/submodules/sdsl-lite/include \
        -I/opt/seqan3/submodules/cereal/include \
        -I/src/Kseq_Cpp \
        main.cpp Kseq_Cpp/Kseq_Cpp.cpp -lz -o /out/IEContig

COPY tests/test_cli.py ./tests/test_cli.py
# A failing regression test fails the image build.
RUN python3 tests/test_cli.py /out/IEContig -v

FROM debian:trixie-slim AS runtime
RUN apt-get update --error-on=any \
    && apt-get install -y --no-install-recommends \
        ca-certificates libgomp1 libstdc++6 zlib1g \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /data \
    && chown 10001:10001 /data

COPY --from=build /out/IEContig /usr/local/bin/IEContig
COPY --from=build /opt/seqan3/LICENSE.md /usr/share/doc/iecontig/seqan3-LICENSE.md
COPY --from=build /opt/seqan3/submodules/sdsl-lite/LICENSE /usr/share/doc/iecontig/sdsl-LICENSE
COPY --from=build /opt/seqan3/submodules/cereal/LICENSE /usr/share/doc/iecontig/cereal-LICENSE
COPY --from=build /opt/seqan3/include/seqan3/contrib/std/LICENSE /usr/share/doc/iecontig/seqan3-std-LICENSE
ENV HOME=/tmp
WORKDIR /data
USER 10001:10001
ENTRYPOINT ["/usr/local/bin/IEContig"]
CMD ["--help"]
