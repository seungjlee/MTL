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
import shutil
import subprocess
import sys
import time

import color
from color import ColorString

# On Windows the default console encoding (cp1252) cannot encode the Unicode
# block characters used by the C++ progress bar. Reconfigure stdout to replace
# unencodable characters instead of crashing the script.
try:
  sys.stdout.reconfigure(errors='replace')
except AttributeError:
  pass

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
parser.add_argument("-Coverage", action='store_true',
                    help='Print test coverage summary after running. Requires gcovr (pip install gcovr) ' +
                    'and a build compiled with --coverage (gcc/clang/MinGW). Not supported with MSVC.')
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

print()
print(ColorString(color.LCYAN, 'Python: ') + sys.version.replace('\n', ''))

if platform.system() == 'Linux':
  import distro
  print(ColorString(color.LCYAN, 'OS: ') + distro.name(pretty=True))
  TestDir = BuildDir + '/Tests/'
else:
  print(ColorString(color.LCYAN, 'OS: ') + platform.platform())
  if args.Debug:
    TestDir = BuildDir + '/Tests/Debug/'
  else:
    TestDir = BuildDir + '/Tests/Release/'

LogFile = BuildDir + '/TestLog_' + datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '.txt'
print(ColorString(color.LCYAN, 'Log File: ') + LogFile)

file = open(LogFile, 'w', encoding="utf-8")

os.chdir(TestDir)
print(ColorString(color.LCYAN, 'Tests Directory: ') + os.getcwd())

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

# Buffer test stdout so we can optionally replay it after the summary,
# but also write each test's output to the log file as soon as it finishes
# so partial logs survive interruptions (Ctrl-C, crash, hang).
Output = io.StringIO()

def WriteToLogAndBuffer(text):
  Output.write(text)
  file.write(text.replace('\r\n', '\n'))
  file.flush()

for test in TestList:
  if (test != '') & (fnmatch.fnmatch(test, Pattern)):
    print('{:.<60}'.format(test), end='')
    sys.stdout.flush()

    testResult = ColorString(color.LGREEN, 'OK ')

    startTime = time.time()
    processResult = subprocess.run(['./' + test] + TestArguments, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, check=False)
    endTime = time.time()

    WriteToLogAndBuffer(TestSeparator.decode('utf-8'))
    WriteToLogAndBuffer(test + '\n')
    WriteToLogAndBuffer(TestSeparator.decode('utf-8'))
    WriteToLogAndBuffer(processResult.stdout.decode('utf-8'))
    WriteToLogAndBuffer('\n\n')

    if processResult.returncode != 0:
      testResult = ColorString(color.RED, 'BAD')
      errorCount = errorCount + 1

    print('{:<7}'.format(testResult), end='')
    print(' %8.3f secs.' % (endTime - startTime,))

WriteToLogAndBuffer(TestSeparator.decode('utf-8'))
totalEndTime = time.time()

if errorCount > 0:
  print(color.RED + '\nNumber of failed tests: %d. Total time: ' % errorCount, end='')
else:
  print(color.LGREEN + '\nAll tests passed. Total time: ', end='')

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

print(color.LCYAN + datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
print(color.RESET)

def PrintCoverageSummary(repoRoot, buildDir):
  gcovr = shutil.which('gcovr')
  if gcovr is None:
    message = "Coverage: gcovr not found. Install with: pip install gcovr"
    print(ColorString(color.LYELLOW, message))
    file.write('\n' + message + '\n')
    file.flush()
    return

  # Tolerate a known gcov bug where heavily-inlined SIMD lines report bogus huge hit
  # counts. See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=68080.
  # --merge-lines collapses per-template-instantiation duplicates so a source line is
  # counted once regardless of how many times the compiler instantiated it. Without
  # this option, line totals balloon with each new template instantiation in the test
  # suite and make the percentage hard to interpret.
  command = [
    gcovr,
    '--root', repoRoot,
    '--filter', os.path.join(repoRoot, 'include') + os.sep,
    '--exclude', os.path.join(repoRoot, 'Tests') + os.sep,
    '--exclude', os.path.join(repoRoot, buildDir) + os.sep,
    # Davis_LDL_COLAMD.h is Timothy Davis's LDL/COLAMD code, vendored in as a
    # third-party header. Excluding it keeps the coverage % focused on
    # first-party MTL code.
    '--exclude', os.path.join(repoRoot, 'include', 'MTL', 'Math', 'Davis_LDL_COLAMD.h'),
    # Test.h is the test framework itself. Its failure-reporting branches
    # (Verify(false), Equal mismatch, exception handlers, ...) only execute
    # when a test actually fails, so they're impossible to cover from a
    # passing test suite. Negative tests live in TestFrameworkNegative.cpp;
    # excluding Test.h here keeps the metric focused on library code.
    '--exclude', os.path.join(repoRoot, 'include', 'MTL', 'Tools', 'Test.h'),
    '--gcov-ignore-parse-errors=suspicious_hits.warn_once_per_file',
    '--merge-lines',
    '--print-summary',
    '--txt',
  ]
  try:
    result = subprocess.run(command, cwd=repoRoot, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, check=False, text=True)
  except OSError as exception:
    message = "Coverage: failed to run gcovr: " + str(exception)
    print(ColorString(color.LYELLOW, message))
    file.write('\n' + message + '\n')
    file.flush()
    return

  if result.returncode != 0:
    message = "Coverage: gcovr exited with code %d." % result.returncode
    print(ColorString(color.LYELLOW, message))
    file.write('\n' + message + '\n')
    if result.stderr:
      print(result.stderr.strip())
      file.write(result.stderr)
    hint = "Build with -DMTL_ENABLE_COVERAGE=ON and re-run the tests so .gcda files are produced."
    print(hint)
    file.write(hint + '\n')
    file.flush()
    return

  # The full --txt report goes to the log file; the trailing --print-summary block (the
  # last 3 lines: lines/functions/branches percentages) goes to stdout.
  print()
  print(ColorString(color.LCYAN, 'Test Coverage:'))
  summaryLines = result.stdout.rstrip().splitlines()[-3:]
  for summaryLine in summaryLines:
    print(summaryLine)

  file.write('\nTest Coverage:\n')
  file.write(result.stdout)
  file.flush()

if args.Coverage:
  PrintCoverageSummary(CurrentDir, BuildDir)

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
    print(line)

sys.exit(errorCount)
