#!/usr/bin/python3

import datetime
import fnmatch
import glob
import math
import os
import subprocess
import sys
import time
import argparse
import platform
import io

totalStartTime = time.time()

parser = argparse.ArgumentParser(description='Run all tests or a subset of tests specified by a pattern.')
parser.add_argument('-p', dest='Pattern', metavar='<pattern>', default='*', type=str,
                    help="Only test names that match this pattern will be run. For example: " +
                    "'Test*'. Default pattern: '*'.")
parser.add_argument('-b', dest='BuildDir', metavar='<build path>',
                    default='Build', type=str,
                    help="Specifies the build directory where the tests are. Default path: 'Build'.")
parser.add_argument("-ConsoleOut", action='store_true', help='Full output to console after summary.')
args = parser.parse_args()

Pattern  = args.Pattern;
BuildDir = args.BuildDir;

SkipTestList = []

TestArguments = ['-NoDisplay', '-DisableProgressBar']
TestSeparatorString = '{:-<100}'.format('')
TestSeparator = TestSeparatorString.encode() + b'\n'
CurrentDir = os.getcwd();

if platform.system() == 'Linux':
  TestDir = BuildDir + '/Tests/'
else:
  TestDir = BuildDir + '/Tests/Release/'

LogFile = BuildDir + '/TestLog_' + datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '.txt'
print('\nLog File: ' + LogFile)

file = open(LogFile, 'w')

os.chdir(TestDir)
print(os.getcwd())

if platform.system() == 'Linux':
  TestList = glob.glob('Test*')
else:
  TestList = glob.glob('Test*.exe')

for skip in SkipTestList:
  for i in range(len(TestList)):
    if fnmatch.fnmatch(TestList[i], skip):
      TestList[i] = ''
      break

TestList.sort()
errorCount = 0
print()

Output = io.StringIO()

for test in TestList:
  if (test != '') & (fnmatch.fnmatch(test, Pattern)):
    print('{:.<60}'.format(test), end='')
    sys.stdout.flush()

    testResult = 'OK'

    startTime = time.time()
    processResult = subprocess.run(['./' + test] + TestArguments, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
    endTime = time.time()

    Output.write(TestSeparator.decode('utf-8'))
    Output.write(processResult.stdout.decode('utf-8'))

    if processResult.returncode != 0:
      testResult = 'FAILED'
      errorCount = errorCount + 1

    print('{:<7}'.format(testResult), end='')
    print(' %8.3f secs.' % (endTime - startTime,))

Output.write(TestSeparator.decode('utf-8'))
totalEndTime = time.time()

if errorCount > 0:
  print('\nNumber of failed tests: %d. Total time: ' % errorCount, end='')
else:
  print('\nAll tests passed. Total time: ', end='')

totalSeconds = totalEndTime - totalStartTime;

if totalSeconds < 60:
  print('%.3f secs.' % totalSeconds)
else:
  totalMinutes = math.floor(totalSeconds / 60);
  seconds = totalSeconds - totalMinutes * 60;
  print('%.0f mins.' % totalMinutes, end='')

  if seconds > 0:
    print(' %.3f secs.' % seconds, end='')

  print(' (%.3f secs.)' % totalSeconds)

if args.ConsoleOut:
  print('')
  print(TestSeparatorString)
  print('Tests Output:')

Output.seek(0)
for line in Output.readlines():
  line = line.replace('\r', '', 1)
  line = line.replace('\n', '', 1)
  file.write(line)
  file.write('\n')
  if args.ConsoleOut:
    print(line)
  
sys.exit(errorCount)
