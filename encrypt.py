import sys
import base64
def enCrypt(message):
    """Totally naive encription routine that will not
    stop a hacker. Use it to encrypt usernames and password.
    Ussage: enCrypt("myusername")"""

    message_bytes = message.encode('ascii')
    base64_bytes = base64.b64encode(message_bytes)
    return base64_bytes.decode('ascii')

if __name__ == "__main__":
   print(sys.argv[1] , "-->", enCrypt(sys.argv[1]))
