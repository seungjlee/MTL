#!/usr/bin/python3
#region pylint disables
# pylint: disable=bad-indentation
# pylint: disable=invalid-name
# pylint: disable=missing-module-docstring
# pylint: disable=consider-using-f-string
# pylint: disable=consider-using-enumerate
# pylint: disable=line-too-long
#endregion

import argparse
import datetime
import fnmatch
import glob
import io
import math
import os
import platform
import subprocess
import sys
import time

import distro

totalStartTime = time.time()

parser = argparse.ArgumentParser(description='Run all tests or a subset of tests specified by a pattern.')
parser.add_argument('-p', dest='Pattern', metavar='<pattern>', default='*', type=str,
                    help="Only test names that match this pattern will be run. For example: " +
                    "'Test*'. Default pattern: '*'.")
parser.add_argument('-b', dest='BuildDir', metavar='<build path>',
                    default='Build', type=str,
                    help="Specifies the build directory where the tests are. Default path: 'Build'.")
parser.add_argument('-Debug', action='store_true', help='Use Debug directory instead of Release directory. (Only on Windows.)')
parser.add_argument("-ConsoleOut", action='store_true', help='Full output to console after summary.')
parser.add_argument("-NoColorRGB24", action='store_true', help='Disable 24-bit RGB colors on test output.')
args = parser.parse_args()

Pattern  = args.Pattern
BuildDir = args.BuildDir

SkipTestList = []

TestArguments = ['-NoDisplay', '-DisableProgressBar']
if args.NoColorRGB24:
  TestArguments.append('-DisableColorRGB24')

TestSeparatorString = '{:-<100}'.format('')
TestSeparator = TestSeparatorString.encode() + b'\n'
CurrentDir = os.getcwd()

print('\033[96m\nOS:\033[0m ' + distro.name(pretty=True))
if platform.system() == 'Linux':
  TestDir = BuildDir + '/Tests/'
else:
  if args.Debug:
    TestDir = BuildDir + '/Tests/Debug/'
  else:
    TestDir = BuildDir + '/Tests/Release/'

LogFile = BuildDir + '/TestLog_' + datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '.txt'
print('\033[96mLog File:\033[0m ' + LogFile)

file = open(LogFile, 'w', encoding="utf-8")

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

    testResult = '\033[92mOK \033[0m'

    startTime = time.time()
    processResult = subprocess.run(['./' + test] + TestArguments, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, check=False)
    endTime = time.time()

    Output.write(TestSeparator.decode('utf-8'))
    Output.write(test + '\n')
    Output.write(TestSeparator.decode('utf-8'))
    Output.write(processResult.stdout.decode('utf-8'))
    Output.write('\n\n')

    if processResult.returncode != 0:
      testResult = '\033[31mBAD\033[0m'
      errorCount = errorCount + 1

    print('{:<7}'.format(testResult), end='')
    print(' %8.3f secs.' % (endTime - startTime,))

Output.write(TestSeparator.decode('utf-8'))
totalEndTime = time.time()

if errorCount > 0:
  print('\033[31m\nNumber of failed tests: %d. Total time: ' % errorCount, end='')
else:
  print('\033[92m\nAll tests passed. Total time: ', end='')

totalSeconds = totalEndTime - totalStartTime

if totalSeconds < 60:
  print('%.3f secs.' % totalSeconds)
else:
  totalMinutes = math.floor(totalSeconds / 60)
  seconds = totalSeconds - totalMinutes * 60
  print('%.0f mins.' % totalMinutes, end='')

  if seconds > 0:
    print(' %.3f secs.' % seconds, end='')

  print(' (%.3f secs.)' % totalSeconds)

print('\033[96m' + datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '\033[0m\n')

if args.ConsoleOut:
  print('')
  print('-------------')
  print('Tests Output:')
  print('-------------')
  print('')

Output.seek(0)
for line in Output.readlines():
  line = line.replace('\r', '', 1)
  line = line.replace('\n', '', 1)
  file.write(line)
  file.write('\n')
  if args.ConsoleOut:
    print(line)
  
sys.exit(errorCount)
