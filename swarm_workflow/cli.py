"""Command line entry point for external swarm workflows."""

from __future__ import annotations

from .commands.parser import build_parser


def main(argv: list[str] | None = None) -> None:
    """Parse one command and delegate it to its owning workflow adapter."""

    parser = build_parser()
    args = parser.parse_args(argv)
    args.handler(args, parser)


if __name__ == "__main__":
    main()
