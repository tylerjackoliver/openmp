#!/bin/bash

# Script to get the project code from the ARCHER system onto the local machine

rsync -azP jtyler@login.archer.ac.uk:/work/y14/y14/jtyler/coursework/affinity_* .
