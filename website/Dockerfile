FROM php:7.4-apache

# Set non-interactive mode for apt-get
ENV DEBIAN_FRONTEND=noninteractive

ARG WEBSITE_INNER_PORT

# Install dependencies for building PHP extensions
RUN apt-get update && apt-get install -y libssl-dev \
    && pecl install mongodb-1.9.1 \
    && docker-php-ext-enable mongodb \
    && apt-get install openssh-client -y \
    && apt-get install -y gridengine-client \
    && apt install r-base -y

# Enable apache modules
RUN a2enmod ssl && a2enmod rewrite

# create .ssh directory for www-data
RUN mkdir -p /var/www/.ssh
RUN chown -R www-data:root /var/www/.ssh

# Configure apache
COPY apache-config.conf /etc/apache2/sites-available/000-default.conf

# PHP ini
COPY custom.ini /usr/local/etc/php/conf.d/custom.ini

# Perform search and replace using sed
RUN sed -i "s/WEBSITE_INNER_PORT/${WEBSITE_INNER_PORT}/g" /etc/apache2/sites-available/000-default.conf

# Change working directory to /app
WORKDIR /app

# Create a cron job file for www-data
RUN echo "0 1 * * * php -f /var/www/html/scripts/clean-cronjob.php >> /logs/cronjobs/cron.log 2>&1" > /etc/cron.d/my-apache-cronjob

# Set the correct permissions for the cron file
RUN chmod 0644 /etc/cron.d/my-apache-cronjob

# Install the cron job for www-data
RUN crontab -u www-data /etc/cron.d/my-apache-cronjob

# Create the entrypoint script: check seed, ssh keys, website code directory, .tmp and cron
RUN echo '#!/bin/sh' > entrypoint.sh && \
    echo 'while [ ! -f /logs/seed/seed-complete ]; do' >> entrypoint.sh && \
    echo '  echo "Waiting for mongo-seed to complete...";' >> entrypoint.sh && \
    echo '  sleep 5;' >> entrypoint.sh && \
    echo 'done;' >> entrypoint.sh && \
    echo 'cp /keys/id_rsa /var/www/.ssh/id_rsa' >> entrypoint.sh && \
    echo 'chmod 600 /var/www/.ssh/id_rsa' >> entrypoint.sh && \
    echo 'chown -R www-data:www-data /var/www/.ssh/id_rsa' >> entrypoint.sh && \
    echo 'set -e' >> entrypoint.sh && \
    echo 'mkdir -p /var/www/html/public' >> entrypoint.sh && \
    echo 'ln -sf /data/Web /var/www/html/public/data' >> entrypoint.sh && \
    echo 'mkdir -p /data/Web/.tmp' >> entrypoint.sh && \
    echo 'chown -R www-data:www-data /data/Web/.tmp' >> entrypoint.sh && \
    echo 'mkdir -p /logs/cronjobs' >> entrypoint.sh && \
    echo 'touch /logs/cronjobs/cron.log' >> entrypoint.sh && \
    echo 'mkdir -p /logs/website' >> entrypoint.sh && \
    echo 'chown -R www-data:www-data /logs/website' >> entrypoint.sh && \
    echo 'cron' >> entrypoint.sh && \
    echo 'exec "$@"' >> entrypoint.sh && \
    chmod +x entrypoint.sh

# Change working directory to /var/www/html
WORKDIR /var/www/html

# Expose the port
EXPOSE ${WEBSITE_INNER_PORT}

# Start the server
ENTRYPOINT ["/app/entrypoint.sh"]
CMD ["apache2-foreground"]