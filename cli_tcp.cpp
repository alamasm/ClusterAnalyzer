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

int main()
{
    //  Create a socket
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
        return 1;
    }

    //  While loop:
    char buf[4096];
    string userInput;

    queue <string> commands;
    do {
        //      Enter lines of text
        cout << "-> ";
        if (commands.empty()) {
            getline(cin, userInput);
            commands.push(userInput);
        }
        else {
            cout << commands.front() << endl;
        }
        userInput = commands.front();
        commands.pop();
        //
        string temp = userInput;
        vector <string> tokens;
        int words = 0;
        size_t pos = 0;
        while ((pos = temp.find(" ")) != string::npos) {
            tokens.push_back(temp.substr(0, pos));
            words++;
            temp.erase(0, pos + 1);
        }
        tokens.push_back(temp);
        words++;
        if (tokens[0] == "run" && words >= 2) {
            ifstream file (tokens[1]);
            while (!file.eof()) {
                getline(file, temp);
                commands.push(temp);
            }
            continue;
        }
        
        
        //      Send to server
        int sendRes = send(sock, userInput.c_str(), userInput.size() + 1, 0);
        if (sendRes == -1)
        {
            cout << "Could not send to server! Whoops!" << endl;
            continue;
        }

        //      Wait for response
        int bytesReceived = recv(sock, buf, 4096, 0);
        if (bytesReceived == -1)
        {
            cout << "There was an error getting response from server" << endl;
        }
        else
        {
            //      Display response
            cout << "SERVER: " << string(buf);
            if (bytesReceived == 1)
                cout << "ok" << endl;
        }
    } while(true);

    //  Close the socket
    close(sock);

    return 0;
}
