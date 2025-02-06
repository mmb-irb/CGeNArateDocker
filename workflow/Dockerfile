FROM alpine:latest

# Install necessary packages
RUN apk add --no-cache bash

# Copy the mock.sh file into the container
COPY mock.sh /usr/local/bin/mock.sh

# Make the mock.sh file executable
RUN chmod +x /usr/local/bin/mock.sh

# Set the entrypoint to execute the mock.sh file
ENTRYPOINT ["/usr/local/bin/mock.sh"]