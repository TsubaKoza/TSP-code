# Simple program to convert numbers to English words up to 999,999

def number_to_words(n):
    if n == 0:
        return "zero"
    if n < 0:
        return "minus " + number_to_words(-n)
    words = []

    units = ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen", "sixteen", "seventeen", "eighteen", "nineteen"]
    tens = ["", "", "twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety"]

    def under_thousand(num):
        result = []
        if num >= 100:
            result.append(units[num // 100])
            result.append("hundred")
            num %= 100
            if num > 0:
                result.append("")
        if num >= 20:
            result.append(tens[num // 10])
            num %= 10
            if num > 0:
                result.append(units[num])
        elif num > 0:
            result.append(units[num])
        return " ".join(filter(None, result))

    if n >= 1000:
        words.append(under_thousand(n // 1000))
        words.append("thousand")
        n %= 1000
        if n == 0:
            return " ".join(words)
        if n < 100:
            words.append("")
    words.append(under_thousand(n))
    return " ".join(filter(None, words))

if __name__ == "__main__":
    num = int(input("Enter a number: "))
    print(number_to_words(num))
