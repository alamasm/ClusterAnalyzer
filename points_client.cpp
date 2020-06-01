#include <iostream>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <unistd.h>
#include <sys/socket.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <string.h>
#include <string>
#include <queue>

using namespace std;

void log(string s, ofstream& log_out) {
    log_out << s << endl;
}

int main() {


    ofstream log_out;
    ifstream in_file;
    log_out.open("out/client_log.txt");
    log("started client", log_out);



        int sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock == -1)
    {
        return 1;
    }

    //  Create a hint structure for the server we're connecting with
    int port = 5555;
    string ipAddress = "127.0.0.1";

    sockaddr_in hint;
    hint.sin_family = AF_INET;
    hint.sin_port = htons(port);
    inet_pton(AF_INET, ipAddress.c_str(), &hint.sin_addr);

    //  Connect to the server on the socket
    int connectRes = connect(sock, (sockaddr*)&hint, sizeof(hint));
    if (connectRes == -1)
    {
        log("error connecting to server", log_out);
        cout << "error connecting to server" << endl;
        return 1;
    }

    //  While loop:
    char buf[4096];
    string userInput;
    
    istream* in;
    int from_file = 1;
    cout << "Read commands from file 1 from konsole 0:" << endl;
    cin >> from_file;
    if (from_file) {
        string filename;
        cout << "Enter filename:" << endl;
        cin >> filename;
        filename = "in/" + filename;
        in_file.open(filename);
        in = &in_file;
    } else {
        cout << "Enter commands here:" << endl;
        log("Reading commands from console", log_out);
        in = &cin;
    }


    

    if (!from_file)getline(*in, userInput);
    while (1) {
        if (!from_file) cout << "command: " << endl;
        getline(*in, userInput);
        //cout << userInput << endl;
        if (userInput == "end") {
            cout << "ending" << endl;
            log("ending", log_out);
            return 0;
        }
        else {
            log("sending request: " + userInput, log_out);
            cout << userInput.c_str() << endl;
            int sendRes = send(sock, userInput.c_str(), userInput.size() + 1, 0);
            if (sendRes == -1)
            {
                cout << "Could not send to server!" << endl;
                log("Could not send to server!", log_out);
                continue;
            }

            //      Wait for response
            cout << "RECIEVING" << endl;
            int bytesReceived = recv(sock, buf, 4096, 0);
            if (bytesReceived == -1)
            {
                cout << "Error getting response" << endl;
                log("Error getting response", log_out);
            }
            else
            {
                //      Display response
                cout << "server: " << string(buf);
                if (bytesReceived == 1)
                    cout << "ok" << endl;
            }
        }
    }
}