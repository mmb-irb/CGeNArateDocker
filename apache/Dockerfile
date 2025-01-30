FROM httpd:2.4

# Copy the custom Apache configuration file
COPY apache-config.conf /usr/local/apache2/conf/conf.d/custom.conf

# Define the build arguments
ARG WEBSITE_INNER_PORT
ARG APACHE_HTTP_INNER_PORT
ARG APACHE_HTTPS_INNER_PORT

# Perform search and replace using sed
RUN sed -i "s/WEBSITE_INNER_PORT/${WEBSITE_INNER_PORT}/g" /usr/local/apache2/conf/conf.d/custom.conf
RUN sed -i "s/APACHE_HTTP_INNER_PORT/${APACHE_HTTP_INNER_PORT}/g" /usr/local/apache2/conf/conf.d/custom.conf
RUN sed -i "s/APACHE_HTTPS_INNER_PORT/${APACHE_HTTPS_INNER_PORT}/g" /usr/local/apache2/conf/conf.d/custom.conf

# Append IncludeOptional directive to the default httpd.conf
RUN echo "IncludeOptional /usr/local/apache2/conf/conf.d/custom.conf" >> /usr/local/apache2/conf/httpd.conf

# Create the SSL directory
# RUN mkdir -p /usr/local/apache2/conf/ssl

# Copy SSL certs into the container if needed
# COPY public.pem /usr/local/apache2/conf/ssl/public.pem
# COPY private.key /usr/local/apache2/conf/ssl/private.key

EXPOSE ${APACHE_HTTP_INNER_PORT} ${APACHE_HTTPS_INNER_PORT}