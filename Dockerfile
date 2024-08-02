# Use a base image that supports multiple architectures
FROM rocker/r-ver:4.1.2

LABEL maintainer="Robert Bowers" \
      name="genome_quality_statistics_app" \
      version="0.1"

# Install system libraries
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libv8-dev \
    pandoc \
    pandoc-citeproc \
    libclang-dev \
    sudo \
    gdebi-core \
    xtail \
    wget \
    gnupg \
    libxt-dev \
    cmake \
    build-essential \
    git \
    r-base-dev

# Install Node.js (required for Shiny Server)
RUN curl -fsSL https://deb.nodesource.com/setup_16.x | bash - && \
    apt-get install -y nodejs

# Install R packages
RUN R -e "install.packages(c('shiny', 'ggplot2', 'broom', 'dplyr', 'DT', 'shinythemes', 'randomcoloR', 'corrr', 'tibble', 'forcats', 'scales'), repos='https://cran.rstudio.com/')"

# Download and install Shiny Server from source
RUN wget https://github.com/rstudio/shiny-server/archive/refs/tags/v1.5.17.973.tar.gz && \
    tar -xzf v1.5.17.973.tar.gz && \
    cd shiny-server-1.5.17.973 && \
    mkdir tmp && cd tmp && \
    cmake .. && make && \
    mkdir -p /usr/local/shiny-server/bin /usr/local/shiny-server/config /usr/local/shiny-server/ext && \
    cp -r ../bin/* /usr/local/shiny-server/bin && \
    cp -r ../config/* /usr/local/shiny-server/config && \
    cp -r ../ext/* /usr/local/shiny-server/ext && \
    cd ../.. && \
    rm -rf shiny-server-1.5.17.973 v1.5.17.973.tar.gz

# Copy necessary files
COPY genome_quality_statistics_app /srv/shiny-server/app

# Make all app files readable (solves issue when running as non-root)
RUN chmod -R 755 /srv/shiny-server/app

# Expose port
EXPOSE 3838

# Start Shiny app
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/app', host = '0.0.0.0', port = 3838)"]
