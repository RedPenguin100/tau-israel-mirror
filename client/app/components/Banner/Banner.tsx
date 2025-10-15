export default function Banner() {
    return (
        <div style={{
            background: "#1a1f36",
            color: "white",
            padding: "10px 14px",
            textAlign: "center",
            fontWeight: 600
        }}>
            This is a demo website  (mirror of the current GitLab) -
            The original website doesn’t work because of one environment variable wasn’t turned on at the wiki freeze. The wiki web interface will be turned on at the Wiki Thaw on November 9.
        </div>
    );
}