FROM ghcr.io/prefix-dev/pixi:0.72.2

LABEL org.opencontainers.image.title="beyondMLST"
LABEL org.opencontainers.image.description="Recombination-aware temporal and contextual bacterial phylogenetics"
LABEL org.opencontainers.image.source="https://github.com/ghruproject/beyondmlst"

WORKDIR /app
COPY . .
RUN pixi install --locked

ENV PYTHONUNBUFFERED=1
ENTRYPOINT ["pixi", "run", "beyondmlst"]
CMD ["--help"]
