
import sys
import socket
import time

HOST = '127.0.0.1'  # Standard loopback interface address (localhost)
PORT = 65432  # Port to listen on (non-privileged ports are > 1023)


def main():
    action = sys.argv[1]
    
    if action == 'start':
        # start socket
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.bind((HOST, PORT))
            sock.listen()
            print(f'Accepting connections under port {PORT}...')


            try:
                while True:
                    connection, address = sock.accept()
                    with connection:
                        operation = connection.recv(1024)

                        if b' ' in operation:
                            (action, param) = operation.split(b' ')
                        else: 
                            action = operation
                            param = b''
                        
                        if action == b'stop':
                            print('Stopping server in 10 seconds...')
                            time.sleep(10)
                            break
                        elif action == b'process':
                            print('Processing BAM...', param)
                        elif action == b'write':
                            print('Write VCF...', param)
                        elif action == b'save':
                            print('Save Checkpoint...')
                        elif action == b'load':
                            print('Load Checkpoint...')
                        

            finally:
                print('Stopping server now...')
                sock.close()




    elif action == 'stop':
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.connect((HOST, PORT))
            sock.sendall(b'stop')
            sock.close()


        print(action)
    elif action == 'process':
        # connect to socket and send add command
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.connect((HOST, PORT))
            sock.sendall(b'process lala.bam')
            sock.close()

    elif action == 'write':
        # connect to socket and send add command
        # connect to socket and send add command
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.connect((HOST, PORT))
            sock.sendall(b'write test.vcf')
            sock.close()
    else:
        print('Invalid Action')
  

if __name__=='__main__':
    main()