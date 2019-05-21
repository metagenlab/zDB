#!/bin/bash
table="$1"
mysql --host=sib-sql04.vital-it.ch -u orthomcl -porthomcl orthomcl -e "drop table $table"
