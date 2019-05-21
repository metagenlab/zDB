#!/bin/bash

mysql --host=sib-sql04.vital-it.ch -u orthomcl -porthomcl orthomcl -e 'SELECT table_schema AS "Database name", SUM(data_length + index_length) / 1024 / 1024 AS "Size (MB)" FROM information_schema.TABLES GROUP BY table_schema;'