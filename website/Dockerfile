FROM php:7.4-apache

# Set non-interactive mode for apt-get
ENV DEBIAN_FRONTEND=noninteractive

ARG WEBSITE_INNER_PORT

# Install dependencies for building PHP extensions
RUN apt-get update && apt-get install -y libssl-dev \
    && pecl install mongodb \
    && docker-php-ext-enable mongodb \
    && apt-get install openssh-client -y \
    && apt-get install -y gridengine-client

# Enable apache modules
RUN a2enmod ssl && a2enmod rewrite
#RUN mkdir -p /etc/apache2/ssl

#COPY ./ssl/*.pem /etc/apache2/ssl/

RUN mkdir -p /root/.ssh
RUN mkdir -p /var/www/.ssh
RUN chown -R www-data:root /var/www/.ssh
# ADD id_rsa /root/.ssh/id_rsa 
# RUN cp /keys/id_rsa /root/.ssh/id_rsa
# RUN chmod 600 /root/.ssh/id_rsa

# Configure apache
COPY apache-config.conf /etc/apache2/sites-available/000-default.conf

# Perform search and replace using sed
RUN sed -i "s/WEBSITE_INNER_PORT/${WEBSITE_INNER_PORT}/g" /etc/apache2/sites-available/000-default.conf

# Change working directory to /app
WORKDIR /app

# Create the entrypoint script to create the symbolic link to the data folder
RUN echo '#!/bin/sh' > entrypoint.sh && \
    # echo 'ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null user@my_stack_sge' >> entrypoint.sh && \
    echo 'cp /keys/id_rsa /root/.ssh' >> entrypoint.sh && \
    echo 'cp /keys/id_rsa /var/www/.ssh' >> entrypoint.sh && \
    echo 'chown -R www-data:www-data /var/www/.ssh/id_rsa' >> entrypoint.sh && \
    echo 'set -e' >> entrypoint.sh && \
    echo 'mkdir -p /var/www/html/public' >> entrypoint.sh && \
    echo 'ln -sf /data/Web /var/www/html/public/data' >> entrypoint.sh && \
    echo 'exec "$@"' >> entrypoint.sh && \
    chmod +x entrypoint.sh

# Change working directory to /var/www/html
WORKDIR /var/www/html

# Expose the port
EXPOSE ${WEBSITE_INNER_PORT}

# Start the server
ENTRYPOINT ["/app/entrypoint.sh"]
CMD ["apache2-foreground"]