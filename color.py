RESET    = "\033[0m" 

BLACK    = "\033[30m"
RED      = "\033[31m"
GREEN    = "\033[32m"
YELLOW   = "\033[33m"
BLUE     = "\033[34m"
MAGENTA  = "\033[35m"
CYAN     = "\033[36m"
WHITE    = "\033[37m"

LBLACK   = "\033[90m"
LRED     = "\033[91m"
LGREEN   = "\033[92m"
LYELLOW  = "\033[93m"
LBLUE    = "\033[94m"
LMAGENTA = "\033[95m"
LCYAN    = "\033[96m"
LWHITE   = "\033[97m"

def ColorString(color, str):
  return color + str + RESET