FROM alpine:latest

# Install necessary packages
RUN apk add --no-cache bash

# Copy the mock.sh file into the container
COPY mock.sh /usr/local/bin/mock.sh

# Make the mock.sh file executable
RUN chmod +x /usr/local/bin/mock.sh

ENV ARG_VALUE=

ENTRYPOINT ["/bin/bash", "-c"]

# Set the entrypoint to execute the mock.sh file
CMD ["/usr/local/bin/mock.sh \"$ARG_VALUE\""]