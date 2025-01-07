import subprocess

class ChessEngine:
    def __init__(self, engine_path):
        self.engine = subprocess.Popen(
            [engine_path],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        self._initialize_engine()

    def _initialize_engine(self):
        self._send_command("uci")
        self._send_command("ucinewgame")
        self._send_command("isready")

    def _send_command(self, command):

        self.engine.stdin.write(command + "\n")
        self.engine.stdin.flush()

    def _read_output(self):
        output = self.engine.stdout.readline().strip()
        return output
    
    def get_best_move(self, fen, movetime, pondering):
        best_move = "nobestmove"
        pondermove = "nopondermove"
        if pondering != "nopondering": 
            self._send_command("stop") #stop pondering
            while True:
                output = self._read_output()
                if output.startswith("bestmove"):
                    parts = output.split()     
                    break
        
        self._send_command(f"position fen {fen}")
        self._send_command(f"go wtime {movetime} btime {movetime} winc 100 binc 100")

        while True:
            output = self._read_output()
            if output.startswith("bestmove"):
                parts = output.split()
                best_move = parts[1]
                pondermove = parts[3] if len(parts) > 3 else "nopondermove"
                break
        
        if movetime < 1500: 
            pondermove = "nopondering"

        if pondermove != "nopondermove" and pondermove != "nopondering":
            self._send_command("go ponder")
        
        return best_move, pondermove

potato2 = None
pondering = "nopondering"

def chess_bot(obs):
    global potato2
    global pondering

    fen = obs['board']
    
    engine_path = '/kaggle_simulations/agent/pototo2'
    #engine_path = '/kaggle/working/pototo2'

    if potato2 is None:
        potato2 = ChessEngine(engine_path)   
        
    time_left = obs['remainingOverageTime'] * 1000  

    best_move, pondering = potato2.get_best_move(fen, time_left, pondering)

    return best_move

