#!/bin/bash

# Script to send the code from the local machine to the ARCHER system

rsync -azP *.f90 /jtyler@login.archer.ac.uk:/work/y14/y14/jtyler/coursework/
