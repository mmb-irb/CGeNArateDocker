FROM php:7.4-apache

ARG WEBSITE_INNER_PORT

# Install dependencies for building PHP extensions
RUN apt-get update && apt-get install -y libssl-dev \
    && pecl install mongodb \
    && docker-php-ext-enable mongodb

# Enable apache modules
RUN a2enmod ssl && a2enmod rewrite
#RUN mkdir -p /etc/apache2/ssl

#COPY ./ssl/*.pem /etc/apache2/ssl/

# Configure apache
COPY apache-config.conf /etc/apache2/sites-available/000-default.conf

# Perform search and replace using sed
RUN sed -i "s/WEBSITE_INNER_PORT/${WEBSITE_INNER_PORT}/g" /etc/apache2/sites-available/000-default.conf

EXPOSE ${WEBSITE_INNER_PORT}

CMD ["apache2-foreground"]