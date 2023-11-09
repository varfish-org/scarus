# Mangement of the GitHub project.

resource "github_repository" "scarus" {
  name        = "scarus"
  description = "Automated evaluation of ACMG rules"

  visibility = "public"

  has_downloads = true
  has_issues = true
  has_projects = false
  has_wiki = false

  allow_rebase_merge = false
  allow_merge_commit = false
}
