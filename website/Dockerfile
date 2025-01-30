# Use nginx Alpine Linux as base image
FROM docker.io/library/nginx:alpine

# Define the build arguments
ARG WEBSITE_INNER_PORT

# Define working dir
WORKDIR /app

# Download and unzip mdposit-client-build repo
RUN git clone https://github.com/mmb-irb/CGeNArateWeb.git

# TODO: edit config.php with env variables and remove config.php from this repo!!!
# RUN cp /app/CGeNArateWeb/app/config.git.php /app/CGeNArateWeb/app/config.php
# RUN echo "DB_HOST=${DB_HOST}" > /app/CGeNArateWeb/app/config.php && \
#     echo "DB_DATABASE=${DB_DATABASE}" >> /app/CGeNArateWeb/app/config.php && \
#     ...
# And then, remove the line below:
COPY config.php /app/CGeNArateWeb/app/config.php

# Copy the built React app to nginx
RUN cp -r /app/CGeNArateWeb/* /usr/share/nginx/html

# Expose port ${WEBSITE_INNER_PORT}
EXPOSE ${WEBSITE_INNER_PORT}

# Start nginx
CMD ["nginx", "-g", "daemon off;"]