#!/bin/sh

while read line ; do
    python master.py --pfamid ${line} --gap_open_cost 50 --force
    sleep 1;
done < family_list
