#!/bin/bash

#Author: Michal Kukowski
#email: michalkukowski10@gmail.com

# This script tests pollard rho

exec=./pollard.out

$exec 2 11 59
$exec 2 424242 5041259
$exec 5 424242 87993167
$exec 7 424242 21441211962585599