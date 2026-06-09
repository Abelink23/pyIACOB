# This file contains the class color, which is use in the
# pyIACOB package to print colored text in the terminal.

class msg:
    # ANSI escape codes for formatting
    COLORS = {
        'p': "\033[95m",   # Purple
        'c': "\033[96m",   # Cyan
        'b': "\033[34m",   # Blue
        'g': "\033[32m",   # Green
        'r': "\033[31m",   # Red
        'y': "\033[33m",   # Yellow
    }

    END = "\033[0m"
    BOLD = "\033[1m"
    UNDER = "\033[4m"
    
    ERROR = "\033[1m" + COLORS['r']  # Bold red for errors
    WARN = "\033[1m" + COLORS['y']   # Bold yellow for warnings
    INFO = COLORS['b']   # Bold blue for info

    def p(self, text: str):
        print(f"{self.COLORS['p']}{text}{self.END}")
    
    def c(self, text: str):
        print(f"{self.COLORS['c']}{text}{self.END}")

    def b(self, text: str):
        print(f"{self.COLORS['b']}{text}{self.END}")

    def g(self, text: str):
        print(f"{self.COLORS['g']}{text}{self.END}")
        
    def r(self, text: str):
        print(f"{self.COLORS['r']}{text}{self.END}")
    
    def y(self, text: str):
        print(f"{self.COLORS['y']}{text}{self.END}")

    def bold(self, color_key: str, text: str):
        """Prints text in bold and the specified color key ('b', 'g', etc.)."""
        color_code = self.COLORS.get(color_key, "")
        print(f"{self.BOLD}{color_code}{text}{self.END}")

    def under(self, color_key: str, text: str):
        """Prints text underlined and the specified color key ('b', 'g', etc.)."""
        color_code = self.COLORS.get(color_key, "")
        print(f"{self.UNDER}{color_code}{text}{self.END}")
    
    def bold_under(self, color_key: str, text: str):
        """Prints text in bold, underlined, and the specified color key ('b', 'g', etc.)."""
        color_code = self.COLORS.get(color_key, "")
        print(f"{self.BOLD}{self.UNDER}{color_code}{text}{self.END}")

    def info(self, text: str):
        print(f"{self.INFO}INFO: {text}{self.END}")
    
    def warn(self, text: str):
        print(f"{self.WARN}WARNING: {text}{self.END}")

    def error(self, text: str):
        print(f"{self.ERROR}ERROR: {text}{self.END}")
