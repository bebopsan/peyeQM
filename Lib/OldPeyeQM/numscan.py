#! /usr/bin/python

## module numscan
''' Una vil copia de la funcion fscanf usada en c
	Permite buscar numeros en una linea de texto.
'''

import re

  

regex = re.compile("(-*\d+\.*\d*|\d*\.\d+)((e|E)(\+|\-)(\d+))*")

  # this object will find numbers within a string

  # usage: regex.findall()

def numscan(fstream):

      '''

      Scan a row of a given filestream

      searching for numbers

      '''

      res = []

      # initialize an empty string

  

      line = fstream.readline()

      # read a line from fstream

  

      nums = regex.findall(line)

      # find numbers in line

  

      for num in nums:

          res.append(eval("%s%s" % (num[0],num[1])))

          # format the number output

          # and append it to the res list

  

      return res 
 
