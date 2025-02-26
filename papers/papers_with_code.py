import os
from dotenv import load_dotenv
from paperswithcode import PapersWithCodeClient
from datetime import datetime

# Load environment variables
load_dotenv()
api_key = os.getenv('PAPERS_WITH_CODE_API_KEY')

# Initialize client with API key
client = PapersWithCodeClient(api_key)

# Get current date for the file name
current_date = datetime.now().strftime('%Y-%m-%d')

# Get papers with correct parameters
papers = client.paper_list(
    q='astrophysics',
    items_per_page=100
)

# Create papers directory if it doesn't exist
os.makedirs('/workspaces/Computational-Astrophysics/papers', exist_ok=True)

# Write to markdown file with date in filename
output_file = f'/workspaces/Computational-Astrophysics/papers/astrophysics_papers_{current_date}.md'

with open(output_file, 'w') as f:
    f.write(f'# Latest Papers on Astrophysics ({current_date})\n\n')
    if papers.results:
        for paper in papers.results:
            title = getattr(paper, 'title', 'No title available')
            authors = ', '.join(getattr(paper, 'authors', [])) or 'No authors available'
            abstract = getattr(paper, 'abstract', 'No abstract available')
            paper_url = getattr(paper, 'paper_url', 'No URL available')
            published = getattr(paper, 'published', 'No date available')

            f.write(f"## {title}\n")
            f.write(f"**Published:** {published}\n")
            f.write(f"**Authors:** {authors}\n")
            f.write(f"**Abstract:** {abstract}\n")
            f.write(f"**URL:** {paper_url}\n\n")
    else:
        f.write('No papers found.')

print(f"Papers have been saved to {output_file}")