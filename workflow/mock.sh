#!/bin/bash

# Set variable from first argument, or default if not provided
ARG_VALUE="${1:-default_value}"

echo "Waiting for the workflow to complete..."

sleep 5

echo "dsfjksdk fjsfhsd fdjksf djksfdhjksf sdf" > /mnt/Web/${ARG_VALUE}/kk

echo "WF completed successfully."
exit 0