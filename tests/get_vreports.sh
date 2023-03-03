#!/usr/bin/env bash

for report in */validation_report.txt; do
  if [ -s $report ]; then
    echo $report
    cat $report
  fi
done
