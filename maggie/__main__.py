import warnings
import typer
import sys
warnings.filterwarnings('ignore')

from .cli import app
#try:
app()
#except Exception as e:
#    typer.echo(f"Error: {e}", err=True)
#    sys.exit(1)