import sys

def query_yes_no(question):
   choices = "[y/n] "
   valid_answer = {'y':True, 'ye':True, 'yes':True, 'n':False, 'ne':False}

   while True:
    sys.stdout.write(question + choices)
    answer = input().lower()
    if answer in valid_answer:
        return valid_answer[answer]
    else:
        sys.stdout.write("Please respond with 'yes' or 'no'. \n")
