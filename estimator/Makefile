init:
	uv sync --all-extras

test:
	uv run pytest tests/ -v

lint:
	uv run ruff check .

format:
	uv run ruff check --fix .
