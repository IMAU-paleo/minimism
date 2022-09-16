#!/usr/bin/env bash

cat << EOF > main.f90
program unit_tests
EOF

for file in test_*.f90; do
  modname=`grep -i module ${file} | grep 'test' | cut -d' ' -f2`
  if [ -n "${modname}" ]; then
    echo "use ${modname}" >> main.f90
  fi
done

for file in test_*.f90; do
  subnames=`grep -i subroutine ${file} | grep 'test' | cut -d' ' -f2`
  for subname in ${subnames}; do
    echo "call ${subname}" >> main.f90
  done
done

cat << EOF >> main.f90
end program unit_tests
EOF
