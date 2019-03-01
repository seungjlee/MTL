import datetime
import fnmatch
import glob
import math
import os
import subprocess
import sys
import time

totalStartTime = time.time()

SkipTestList = []

TestArgument0 = '-NoDisplay'
TestSeparator = '{:-<80}\n'.format('').encode()
CurrentDir = os.getcwd();
TestDir = 'Build/Tests/Release/'
Pattern = '*'

LogFile = CurrentDir + '\\TestLog_' + datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '.txt'
print('\nLog File: ' + LogFile)

file = open(LogFile, 'wb')

NumberOfArguments = len(sys.argv)

if NumberOfArguments == 2:
  Pattern = sys.argv[1]

os.chdir(TestDir)

TestList = glob.glob('Test*.exe')

for skip in SkipTestList:
  for i in range(len(TestList)):
    if fnmatch.fnmatch(TestList[i], skip):
      TestList[i] = ''
      break

errorCount = 0
print()
for test in TestList:
  if (test != '') & (fnmatch.fnmatch(test, Pattern)):
    print('{:.<60}'.format(test), end='')
    sys.stdout.flush()

    testResult = 'OK'

    startTime = time.time()
    # Python 3.3
    try:
      processResult = subprocess.check_output([test, TestArgument0], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
      processResult = e.output
      testResult = 'FAILED'
      errorCount = errorCount + 1
    # Python 3.5
    # processResult = subprocess.run(test + TestArguments, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
    endTime = time.time()

    file.write(TestSeparator)
    # Python 3.3
    file.write(processResult)
    # Python 3.5
    # file.write(processResult.stdout)
    # if processResult.returncode != 0:
    #   testResult = 'FAILED'
    #   errorCount = errorCount + 1

    print('{:<7}'.format(testResult), end='')
    print(' %8.3f secs.' % (endTime - startTime,))

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

sys.exit(errorCount)
