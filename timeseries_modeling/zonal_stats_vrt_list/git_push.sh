# A simple script to add, commit, and push
# Usage: ./git_push.sh "your commit message"

if [ -z "$1" ]; then
  echo "Usage: $0 'commit message'"
  exit 1
fi

git add .
git commit -m "$1"
git push