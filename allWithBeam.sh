#!/bin/bash
for run in 154 155 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176
do
  echo "#####  316 " $run "   #####"
  source dqmLxbatch.sh 316 $run
done
