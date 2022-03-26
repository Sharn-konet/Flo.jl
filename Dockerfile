FROM julia:1.7.2-alpine AS base

LABEL MAINTAINER="Sharn-Konet Reitsma"

WORKDIR /workdir/

COPY . .

RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'


FROM base as application

CMD ["julia", "App.jl"]