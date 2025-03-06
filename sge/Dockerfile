# Use an official Debian image as the base image
FROM debian:bullseye-20241111-slim

# Set non-interactive mode for apt-get
ENV DEBIAN_FRONTEND=noninteractive

ARG DOCKER_GROUP_ID
ARG SUBMITTER_HOSTNAME
ARG SGE_SSH_PORT
ARG SGE_PORT
ARG SGE_MAX_JOBS

# Install necessary packages: gnupg, vim, parallel, g++
RUN apt-get update -y  \
    && apt-get install -y \
        gnupg2 \
        vim \
        parallel \
        g++ \
        docker.io \
        openssh-server \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy the local keys into the Docker image
COPY keys /keys

# Import GPG keys
RUN gpg --import /keys/AA8E81B4331F7F50.asc \
    && gpg --import /keys/DCC9EFBF77E11517.asc \
    && gpg --import /keys/112695A0E562B32A.asc \
    && gpg --import /keys/04EE7237B7D453EC.asc \
    && gpg --import /keys/EF0F382A1A7B6500.asc \
    && gpg --import /keys/648ACFD622F3D138.asc \
    && gpg --import /keys/DCC9EFBF77E11517.asc \
    && gpg --import /keys/54404762BBB6E853.asc \
    && gpg --import /keys/0E98404D386FA1D9.asc

# Update the package lists and install GridEngine packages
RUN apt-get update -y && apt-get install -y gridengine-exec gridengine-master gridengine-qmon

# Copy the public key to the authorized_keys file
RUN mkdir -p /var/www/.ssh
RUN chown -R www-data:www-data /var/www

# Make sure host and container share the same GID for group 'docker', bc it has reading permissions to the socket file
RUN groupmod -g $DOCKER_GROUP_ID docker
RUN usermod -s /bin/bash www-data
RUN usermod -aG docker www-data

# Add setup script and set permissions
ADD setup_gridengine.sh /usr/local/bin/setup_gridengine.sh
ADD setup_submitter.sh  /usr/local/bin/setup_submitter.sh
RUN chmod +x /usr/local/bin/setup_gridengine.sh
RUN chmod +x /usr/local/bin/setup_submitter.sh

# Add the submitter hostname to the submit_hosts.list file
RUN echo "$SUBMITTER_HOSTNAME" > /var/www/submit_hosts.list

# Change the SSH port from 22 default to the one specified in the environment variable
ENV SGE_SSH_PORT=${SGE_SSH_PORT}
RUN sed -i "s/#Port 22/Port ${SGE_SSH_PORT}/" /etc/ssh/sshd_config

# Expose SSH and SGE ports
EXPOSE ${SGE_SSH_PORT} ${SGE_PORT}

# Set the maximum number of jobs as an environment variable for making it accessible in setup_gridengine.sh
ENV SGE_MAX_JOBS=${SGE_MAX_JOBS}

WORKDIR /app

# Create the entrypoint script: check seed, ssh keys and /var/run/docker.sock
RUN echo '#!/bin/sh' > entrypoint.sh && \
    echo 'while [ ! -f /logs/seed/seed-complete ]; do' >> entrypoint.sh && \
    echo '  echo "Waiting for mongo-seed to complete...";' >> entrypoint.sh && \
    echo '  sleep 5;' >> entrypoint.sh && \
    echo 'done;' >> entrypoint.sh && \
    echo 'cp /keys/id_rsa.pub /var/www/.ssh/authorized_keys' >> entrypoint.sh && \
    echo 'chown www-data:www-data /var/www/.ssh/authorized_keys' >> entrypoint.sh && \
    echo 'service ssh start' >> entrypoint.sh && \
    echo 'chmod 666 /var/run/docker.sock' >> entrypoint.sh && \
    echo 'exec "$@"' >> entrypoint.sh && \
    chmod +x entrypoint.sh

ENTRYPOINT ["/app/entrypoint.sh"]
CMD ["bash", "-c", "/usr/local/bin/setup_gridengine.sh && /usr/local/bin/setup_submitter.sh /var/www/submit_hosts.list && qconf -mconf && tail -f /var/spool/gridengine/qmaster/messages"]

# CMD ["bash", "-c", "cp /keys/id_rsa.pub /var/www/.ssh/authorized_keys && chown www-data:www-data /var/www/.ssh/authorized_keys && service ssh start && chmod 666 /var/run/docker.sock && cat /var/www/submit_hosts.list && /usr/local/bin/setup_gridengine.sh && cat /var/www/submit_hosts.list && /usr/local/bin/setup_submitter.sh /var/www/submit_hosts.list && qconf -mconf && tail -f /var/spool/gridengine/qmaster/messages"]